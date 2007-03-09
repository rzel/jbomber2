package game;

import java.net.*;
import java.io.*;

import network.*;
import config.*;
import network.*;

/**
 * A <code>DelegatorClient</code> recieves commands from a <code>
 * NetworkClient</code> and handles this data in a proper way. A <code>
 * DelegatorClient</code> also recieves data from <code>PlayerClient</code>,
 * This data is send to the connected server.
 */
public class DelegatorClient extends Thread
{
    // NetworkClient to communicate with server
    private NetworkClient network_client;

    // The GameEngine that is used to play this game
    private GameEngine game_engine;

    // The id of the player that uses this DelegatorClient
    // -1 indicates that the server has not yet send a player id.
    private int player_id = -1;

    // The nickname of this player
    private String nickname;

    // The command that needs to be send each interval
    private Command to_send;

    /**
     * Creates a new <code>DelegatorClient</code> object, using a given
     * <code>NetworkCommand</code> to communicate with the connected
     * server.
     *
     * @param network_client The NetworkClient that is used to communicate
     *                       with the server.
     * @param nickname       The nickname of this player.
     */
    public DelegatorClient(InetAddress s_addr, int s_port,
                           PrintStream log, String nickname)
    {
        this.network_client = new NetworkClient(s_addr, s_port, this, log);
        this.nickname       = nickname;
        start();
    }

    /**
     * Sends <code>Command</code>s every interval step, if no
     * <code>Command</code> was specified a <code>PlayerCommand</code>
     * will be sent. The <code>getActionID()</code> method of this
     * command will be <code>Config.ACT_NOTHING</code>.
     */
    public void run()
    {
        for (;;)
        {
            // TODO vanalles

            // Sleep for 20 ms.
            try
            {
                sleep(50);
            }
            catch (InterruptedException ie) {}

            // if no command was specified during this sleep,
            // send the NOTHING PlayerCommand, else send the
            // specified command
            if (to_send == null)
            {
                to_send = new PlayerCommand(player_id, Config.ACT_NOTHING);
            }

            network_client.send(to_send);
            to_send = null;
        }
    }

    /**
     * Returns whether it is possible for this client to handle more data.
     * It is possible to handle data, when no other data is to be send during
     * the current interval step.
     *
     * @return <code>true</code> if data can be handled, <code>false</code>
     *         otherwise.
     */
    public boolean canHandleData()
    {
        return to_send == null;
    }

    /**
     * Handles an action of the player that uses this client.
     *
     * @param action_id The id of the action to handle.
     */
    public void handleAction(int action_id)
    {
        // Only handle actions when it is possible
        if (!canHandleData()) return;

        // Make sure only known actions are handled
        switch (action_id)
        {
            case Config.ACT_MOVE_RIGHT:
            case Config.ACT_MOVE_LEFT:
            case Config.ACT_MOVE_DOWN:
            case Config.ACT_MOVE_UP:
            case Config.ACT_LEAVE:
                to_send = new PlayerCommand(player_id, action_id);
        }
    }

    /**
     * Hanldes incomming <code>Command</code>s from the
     * <code>NetworkCommand</code> that is connected to a server.
     *
     * @param command The command to handle.
     */
    public void handleCommand(Command command)
    {
        // Handle all distinct commands properly.
        // Casting should not be a problem.
        switch (command.getCommandID())
        {
            case Config.CMD_GAME:
                // Handle GameCommand
                GameCommand gcommand = (GameCommand)command;

                // Handle all the commands in a GameCommand
                for (Command c : gcommand.getCommands())
                {
                    handleCommand(c);
                }

                break;
            case Config.CMD_INIT:
                // Handle an InitCommand.
                InitCommand icommand = (InitCommand)command;
                // Make the link a player and this DelegatorClient
                player_id = icommand.getClientID();

                // Initialize game engine
                if (game_engine == null)
                {
                    // Try to initialize game_engine
                    try
                    {
                        game_engine = new GameEngine(icommand.getSchemeName());
                    }
                    catch (Exception e)
                    {
                        // Initializing game_engine failed
                        System.out.println("Error while loading scheme.");
                        System.exit(-1);
                    }
                }
                else
                {
                    // TESTING
                    System.out.println("GameEngine already initialized.");
                    System.exit(-1);
                }

                // Create NewPlayerCommand with true flag to send to server
                // with true flag, so everyone knows this player is new
                int x_pos = icommand.getXPos();
                int y_pos = icommand.getYPos();
                NewPlayerCommand npc = new NewPlayerCommand(player_id,
                                                            x_pos,
                                                            y_pos,
                                                            nickname,
                                                            true);
                // Send NewPlayerCommand to server
                to_send = npc;
                System.out.println("" + npc + " sent to server.");
                break;
            case Config.CMD_NEW:
                // Handle a NewPlayerCommand
                NewPlayerCommand npcommand = (NewPlayerCommand)command;

                // Add player to game_engine when it is not yet in the game
                if (game_engine.getPlayer(npcommand.getSenderID()) == null)
                {
                    game_engine.addPlayer(npcommand.getSenderID(),
                                          npcommand.getNickname(),
                                          npcommand.getXPos(),
                                          npcommand.getYPos());
                }

                // If this was a new player, send command to notify this player
                // that this client exists.
                if (npcommand.getNewPlayer())
                {
                    // The player that plays using this client
                    ReadablePlayer p = game_engine.getPlayer(player_id);
                    npc = new NewPlayerCommand(player_id,
                                               p.getXPos(),
                                               p.getYPos(),
                                               p.getNickname(),
                                               false);

                    // Send NewPlayerCommand to server
                    to_send = npc;
                }

                break;
            case Config.CMD_PLAYER:
                // Handle PlayerCommand
                PlayerCommand pcommand = (PlayerCommand)command;

                game_engine.handleAction(pcommand.getSenderID(),
                                         pcommand.getActionID());
                break;
        }
    }

    /**
     * Closes the connection to the server.
     */
    public void close()
    {
        network_client.send(new PlayerCommand(player_id, Config.ACT_LEAVE));
        network_client.close();
    }

    /**
     * Returns all the <code>ReadablePlayers</code> that play this game.
     *
     * @return all the players that play this game.
     */
    public ReadablePlayer[] getPlayers()
    {
        return game_engine.getPlayers();
    }

    /**
     * Returns the scheme that is used in this game.
     *
     * @return The scheme that is used by this game, if it's spceified,
     *         else <code>null</code>.
     */
    public char[][] getScheme()
    {
        return game_engine.getScheme();
    }

}