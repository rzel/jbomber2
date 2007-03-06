import java.io.Serializable;

/**
 * A <code>Command</code> is send by a client to a server,
 * the server forwards it in a way to all other clients.
 * Every <code>Command</code> has an unique id, so commands
 * can be disinguished by means of this id.
 * A <code>Command</code> also contains a "sender id", so
 * is possible to know what client send the command.
 * Besides these two id's a <code>Command</code> can hold
 * command specific information.
 */
public class Command implements Serializable
{
    // Something
    public static final long serialVersionUID = 1;

    // The unique id to distinguish this command.
    private int command_id;

    // The unique id to distinguish the sender of the command.
    private int sender_id;

    /**
     * Creates a new <code>Command</code> object, whith the given
     * command id and the given sender id.
     *
     * @param command_id The unique id to distinguish this command.
     * @param sender_id  The unique id to distinguish the sender.
     */
    public Command(int command_id, int sender_id)
    {
        this.command_id = command_id;
        this.sender_id  = sender_id;
    }

    /**
     * Returns the unique id to distinguish this command.
     *
     * @return The unique id to distinguish this command.
     */
    public int getCommandID()
    {
        return this.command_id;
    }

    /**
     * Returns the unique id to distinguish the sender.
     *
     * @return The unique id to distinguish the sender.
     */
    public int getSenderID()
    {
        return this.sender_id;
    }
}