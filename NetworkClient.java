import java.io.*;
import java.net.*;

/**
 * A <code>NetworkClient</code> connects is connected to a <code>Server</code>.
 * Commands are send between the server and the client.
 * A <code>NetworkClient</code> doesn't handle incomming commands, it only
 * servers a sender/reciever.
 */
public class NetworkClient extends Thread
{
    // The socket that is connected to a server
    private Socket server;

    // Stream to write log messages to
    private PrintStream log;

    // ObjectOutputStream to send commands
    private ObjectOutputStream obj_out;

    // ObjectInputStream to recieve commands
    private ObjectInputStream obj_in;

    // Indicates whether this server is running.
    private boolean running = true;

    // The DelegatorClient that uses this NetworClient
    private DelegatorClient delegator_client;

    /**
     * Creates a new <code>NetworkClient</code> object.
     *
     * @param s_addr The <code>InetAddress</code> of a server.
     * @param s_port The port on which a server is running.
     * @param user   The <code>DelegatorClient</code> that uses this object.
     * @param log    The <code>PrintStream</code> to write log messages to.
     */
    public NetworkClient(InetAddress s_addr, int s_port,
                         DelegatorClient delegator_client, PrintStream log)
    {
        this.log = log;
        this.delegator_client = delegator_client;

        // Try to connect and initialize in and output streams
        try
        {
            // Try to connect
            server = new Socket(s_addr, s_port);

            // Try to initialize obj_in
            InputStream in = server.getInputStream();
            obj_in = new ObjectInputStream(in);

            // Try to initialize obj_out
            OutputStream out = server.getOutputStream();
            obj_out = new ObjectOutputStream(out);

            writeLog("Client connected!");
            start();
        }
        catch (IOException ioe)
        {
            // Error while connection or initializing.
            writeLog("Error while connecting, message: " + ioe.getMessage());
            // Close connection
            close();
        }
    }

    /**
     * Recieves <code>Commands</code> from a connected server.
     * If something goes wrong, the connection is closed.
     */
    public void run()
    {
        while (running)
        {
            // Try to recieve commands
            try
            {
                Command command = (Command)obj_in.readObject();
                delegator_client.handleCommand(command);
            }
            catch (IOException ioe)
            {
                // Thrown if client was blocking on closing
                writeLog("Error while recieving commands, message " +
                         ioe.getMessage());
                close();
            }
            catch (ClassNotFoundException cnfe)
            {
                // Thrown if client was blocking on closing
                writeLog("Error while recieving commands, message " +
                         cnfe.getMessage());
                // Thrown when recieving data failed.
                // Close this connection.
                close();
            }

            // Sleep for a while
            try
            {
                sleep(20);
            }
            catch (InterruptedException ie) {}
        }
    }

    /**
     * Sends data to a connected server. If sending fails,
     * this <code>NetworkClient</code> is closed.
     */
    public void send(Command command)
    {
        try
        {
            if (running) obj_out.writeObject(command);
        }
        catch (IOException ioe)
        {
            writeLog("Error while sending, message" + ioe.getMessage());
            close();
        }
    }

    /**
     * Returns the closed state of the </code>NetworkClient</code>.
     *
     * @Return <code>true</code> if the <code>NetworkClient</code> has been closed.
     */
    public boolean isClosed()
    {
        return server.isClosed();
    }

    /**
     * Closes the connection to the connected server.
     * Sending and recieving commands becomes impossible.
     */
    public void close()
    {
        running = false;

        // Try to close connection
        try
        {
            if (server != null)
            {
                server.close();
            }
        }
        catch (IOException ioe)
        {
             // Thrown if client was blocking on closing
            writeLog("Error while closing connection, message " +
                     ioe.getMessage());
        }
    }

    /**
     * Writes log messages to log if it exsits.
     *
     * @param message The message to write to log.
     */
    private void writeLog(String message)
    {
        if (log != null)
        {
            log.println(message);
        }
    }
}